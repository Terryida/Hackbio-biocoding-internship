details <- list(
  first_name = c("Tyler", "Mahpara", "Peter", "Gracious", "Terry"),
  last_name = c("Maire", "Abid", "Imonte", "Isoah", "Ngala"),
  slack_username = c("Tyler", "Mahpara", "Perex", "Gracious", "Terry"),
  email = c("tam357a@gmail.com", "mahparaabid97@gmail.com", "imontepez@gmail.com", "isoahgracious@gmail.com", "terryngala9@gmail.com"),
  hobby = c("Carpentry", "Enjoying nature", "Playing music", "Reading", "Travelling"),
  country = c("United States", "Pakistan", "Nigeria", "Nigeria", "Kenya"),
  discipline = c("Bioinformatics", "Biotechnology", "Biomedical Technology", "Microbiology", "Biochemistry"),
  programming_language = c("R", "R", "R", "R", "R")
)
#insert the list above into a data frame then print
Interns_data <- data.frame(
  FirstName = details$first_name,
  LastName = details$last_name,
  Slack_UserName = details$slack_username,
  Email = details$email,
  Hubby = details$hubby,
  Country = details$country,
  Discipline = details$discipline,
  Programming_Language = details$programming_language
)
print(Interns_data)
# Group Github Usernames
#https://github.com/mahpara97, https://github.com/Terryida
#https://github.com/imontepez , https://github.com/tylermaire , https://github.com/Graciousisoah
